Parallel Matrix Multiplication and Inversion using MPI
This project explores the implementation and performance analysis of parallel matrix multiplication and matrix inversion using the Message Passing Interface (MPI) standard. It supports both dense and sparse matrices, with a focus on scalability, communication efficiency, and computational performance in distributed environments.

📌 Key Features
✅ Parallel dense matrix multiplication using MPI

✅ Parallel sparse matrix multiplication using CSR (Compressed Sparse Row) format

✅ Parallel matrix inversion for dense matrices using Gauss-Jordan elimination

📊 Performance evaluation under strong and weak scalability scenarios

☁️ Experiments conducted across various cluster types (Fat, Light, Infra-regional)

🔧 Technologies Used
C++

MPI (OpenMPI / MPICH)

Google Cloud Platform for VM-based cluster simulation

CSR (Compressed Sparse Row) format for sparse matrix efficiency

📈 Performance Metrics
The project evaluates and compares the following:

Execution time

Speedup (compared to serial implementation)

Communication overhead

Efficiency and scalability (using Amdahl’s Law)

📊 Experimental Setup
Experiments were conducted on three cluster types:

Fat Cluster (intra-regional): High-performance VMs for low-latency testing

Light Cluster: Resource-constrained VMs to test behavior under limited resources

Infra-regional Cluster: Distributed VMs across regions to evaluate communication latency

Sparse matrix operations used a fixed sparsity of 10% for consistent benchmarking.

⚠️ Limitations
Sparse matrix inversion was not implemented directly due to its high complexity and low scalability.

Inversion operations show less speedup due to their sequential nature and frequent synchronization requirements.

📚 Report
A complete academic report with theoretical analysis, performance graphs, and conclusions is included in the repository under /Report/.

👩‍💻 Author
Yasamin Hosseinzadeh Sani
Department of Computer Engineering - Data Science
University of Pavia, Italy
Contact: yasamin.hosseinzadehsa01@universitadipavia.it
