#include <iostream>

#include <omp.h>
#include "vector"
#include "cmath"


void printArray(std::vector<int> arr){
    for (int j=0; j < arr.size(); ++j) {
        std::cout << arr[j] << " ";
    }
    std::cout << "\n";
}


void unScalablePrefixSum(std::vector<int>& input){
    // Number of threads equals the size of input array
    int size = input.size();
    omp_set_num_threads(size);
    std::vector<int> helpRegister(size);

    // Output array initialization
    # pragma omp parallel for shared(size, helpRegister, input) default(none)
    for (int i = 0; i < size; ++i) {
        printf("hello from unScalablePrefixSum from thread %d\n", omp_get_thread_num());
        helpRegister[i] = input[i];
    }


    for (int j = 0; j < log2(size); ++j) {
        # pragma omp parallel for shared(helpRegister, input, size, j) default(none)
        for (int k = pow(2, j); k < size; ++k) {
            helpRegister[k] = helpRegister[k] + input[k - pow(2, j)];
        }
        # pragma omp parallel for shared(helpRegister, input, size, j) default(none)
        for (int h = pow(2, j); h < size; ++h) {
            input[h] = helpRegister[h];
        }
    }
}


std::vector<int> scalablePrefixSum(std::vector<int> input, int num_threads){
    // variable initialization
    std::vector<int> zVector(num_threads);
    std::vector<int> sVector(input.size());

    int chunk_size = input.size()/num_threads;
    int k;
    int sum;
    int index;

    omp_set_num_threads(num_threads);
    # pragma omp parallel for shared(num_threads, zVector, sVector, chunk_size, input) private(sum, k, index) default(none)
    for (int i = 0; i < num_threads; ++i) {
        sum = 0;
        printf("hello scalablePrefixSum from thread %d\n", omp_get_thread_num());

        // sequential computing of prefix sum for every thread sub-vector
        for (k = 0; k < chunk_size; ++k) {
            index = (i*chunk_size) + k;
            sum = sum + input[index];
            sVector.at(index) = sum;
        }
        zVector.at(i) = sum;
    }

    // Call unscalable prefix sum for zVector variable
    unScalablePrefixSum(zVector);

    omp_set_num_threads(num_threads);
    # pragma omp parallel for shared(zVector, num_threads, chunk_size, sVector) private(k, index) default(none)
    for(int i = 1; i < num_threads; ++i){
        for (k=0; k < chunk_size; ++k){
            index = (i*chunk_size) + k;
            sVector[index] = sVector[index] + zVector[i-1];
        }
    }

    return sVector;

}


int main() {

    int algorithmChoice;
    std::cout << "Select the algorithm:\n";
    std::cout << "1. EREW PRAM algorithm\n";
    std::cout << "2. Scalable EREW PRAM algorithm\n";
    std::cin >> algorithmChoice;

    std::vector<int> output;


    if (algorithmChoice == 1) {

        int arraySize;
        std::cout << "Enter the array size (Max size is 32): ";
        std::cin >> arraySize;
        std::vector<int> input(arraySize);

        // Check array size
        if (arraySize > 32) {
            std::cout << "Array size exceeded the limit (32) for EREW PRAM algorithm." << std::endl;
            return 0;
        }

        // Automatically generate input array values as a sequence of 1 - min(N, 32)
        for (int i = 0; i < arraySize; i++) {
            input[i] = i + 1;
        }

        std::cout << "This is generated input array:\n";
        printArray(input);

        unScalablePrefixSum(input); // Calling alg.
        std::cout << "Result of unscalable EREW PRAM prefix sum:\n";
        printArray(input);
        return 1;

    } else if(algorithmChoice == 2){

        int arraySize;
        std::cout << "Enter the array size: ";
        std::cin >> arraySize;
        std::vector<int> input(arraySize);

        int numThreads;
        std::cout << "Enter the number of threads: ";
        std::cin >> numThreads;

        if (arraySize % numThreads != 0){
            std::cout << "Please note that if array size modulo number of threads does not equal to 0, the array cannot be split "
                         "into equal chunks, which may result in unexpected behavior of this algorithm.\n";
        }

        // Automatically generate input array values as a sequence of 1 - min(N, 32)
        for (int i = 0; i < arraySize; i++) {
            input[i] = i + 1;
        }

        printArray(input);
        std::vector<int> res = scalablePrefixSum(input, numThreads);
        printArray(res);
        return 1;
    } else {
        std::cout << "Invalid input";
        return 1;
    }
}
