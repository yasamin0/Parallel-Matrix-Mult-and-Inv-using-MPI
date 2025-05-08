# 🧮 Parallel Matrix Multiplication and Inversion using MPI

This project explores the implementation and performance analysis of **parallel matrix multiplication** and **matrix inversion** using the [Message Passing Interface (MPI)](https://www.mpi-forum.org/).  
It supports **dense** and **sparse** matrices, focusing on scalability, communication efficiency, and computational performance across distributed environments.

---

## 🔑 Key Features

- ✅ **Parallel dense matrix multiplication** using MPI
- ✅ **Parallel sparse matrix multiplication** using CSR (Compressed Sparse Row) format
- ✅ **Parallel matrix inversion** (for dense matrices) using Gauss-Jordan elimination
- 📊 Performance evaluation under **strong and weak scalability** scenarios
- ☁️ Experiments on various cluster types: **Fat**, **Light**, and **Infra-regional**

---

## 🧰 Technologies Used

- C++  
- MPI (OpenMPI / MPICH)  
- CSR format for sparse matrices  
- Google Cloud Platform (GCP) for cluster deployment  

---

## 📁 Project Structure

/Code/
├── serial_mult.cpp # Serial matrix multiplication (dense & sparse)
├── serial_inv.cpp # Serial matrix inversion (Gauss-Jordan)
├── parallel_mult.cpp # Parallel matrix multiplication with MPI
├── parallel_inv.cpp # Parallel matrix inversion with MPI
└── /scripts/ # Shell scripts for automation and testing


---

## ⚙️ Cluster Configurations

| Cluster Type   | VM Specs                    | Location         |
|----------------|-----------------------------|------------------|
| **Fat**        | High-performance (n2-standard-16) | Intra/Infra-regional |
| **Light**      | Low-resource (e2-medium/small)   | Intra-regional     |
| **Infra**      | Distributed across zones     | Inter-regional    |

---

## 📈 Performance Metrics

- Execution time
- Speedup (vs. serial)
- Communication overhead
- Efficiency (using Amdahl’s Law)

---

## 📄 Report

For in-depth theory, methodology, and results, see the full report in the [`/Report/`](./Report) folder.  
Includes graphs, scalability analysis, and cluster performance comparisons.

---

## 👩‍💻 Author

**Yasamin Hosseinzadeh Sani**  
Department of Computer Engineering - Data Science  
University of Pavia, Italy  
📧 yasamin.hosseinzadehsa01@universitadipavia.it

---

## 📜 License

This project is licensed under the **MIT License**. See the [LICENSE](./LICENSE) file for details.
