#!/usr/bin/env python3
import subprocess

tests = ["./test_02karatsuba", "./test_03toom", "./test_04ntt", "./test_05fft",
         "./test_06radix3fft", "./test_07incomplete", "./test_08goods"]
iterations = 10

for test in tests:
    for _ in range(iterations):
        subprocess.check_output([test])
    print(f"ran {iterations} iterations of {test} test. ALL GOOD.")