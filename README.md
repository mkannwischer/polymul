This code accompanies my PhD thesis on **Polynomial Multiplication for Post-Quantum Cryptography**.

If you find mistakes, parts that could use some more comments, examples that could be added, or anything else that you think should be improved, please drop me an e-mail at `<matthias@kannwischer.eu>`. 
I'm also planning to extend this collection to more algorithms that others have used to implement PQC. If you have input about what should be covered, I'd be grateful for any input.

**What this code tries to achieve**

- Demonstrate the core ideas of polynomial multiplication implementations over the last years in a readable and easy to understand fashion.

**What this code does not try to achieve**

- Setting new speed records at anything. I tried to keep the code general and whenever possible independent of the used parameters. In particular, I tried to keep orthogonal optimization strategies separate to make it easier to grasp them. This is of course not yielding an optimal implementation.
- Provide any security guarantees at all. A lot of the code is not implemented in constant time (especially those parts that are usually pre-computed anyway). Hence, this code is not production-ready and before using it one needs to eliminate all the timing leaks. Additionally, I did not extensively verify correctness for a wide range of parameters. There are likely parameters that seem to be supported but cause some overflows in some places which in the best case renders them completely non-functional and in the worst-case works most of the time but fails for extreme values. Careful range analysis is required for actual parameters.



# Structure

| Thesis Chapter       | Python Example                                | C Example                             | 
| -------------------- | --------------------------------------------- | ------------------------------------- | 
| 2.2.1 Schoolbook     | [01schoolbook.py](./Python/01schoolbook.py)   | [01schoolbook.c](./C/01schoolbook.c)  |
| 2.2.2 Karatsuba      | [02karatsuba.py](./Python/02karatsuba.py)     | [02karatsuba.c](./C/02karatsuba.c)    |
| 2.2.3 Toom--Cook     | [03toom.py](./Python/03toom.py)               | [03toom.c](./C/03toom.c)              |
| 2.2.4 NTT            | [04ntt.py](./Python/04ntt.py)                 | [04ntt.c](./C/04ntt.c)                |
| 2.2.5 Radix-2 FFT    | [05fft.py](./Python/05fft.py)                 | [05fft.c](./C/05fft.c)                |
| 2.2.6 Radix-3 FFT    | [06radix3fft.py](./Python/06radix3fft.py)     | [06radix3fft.c](./C/06radix3fft.c)    |
| 2.2.7 Incomplete NTT | [07incomplete.py](./Python/07incomplete.py)   | [07incomplete.c](./C/07incomplete.c)  |
| 2.2.8 Good's trick   | [08goods.py](./Python/08goods.py)          | [08goods.c](./C/08goods.c)          |

Note that the NTT implementations do modular reductions using `% q`. A real implementation would not do that since it is both slow and potentially vulnerable to timing attacks. In a real implementation one would use Montgomery or Barrett reductions. This also affects the twiddle factors as those need to be transformed to Montgomery domain before. 

For more details on the general framework, see [Python/README.md](./Python/README.md) and [C/README.md](./C/README.md]

