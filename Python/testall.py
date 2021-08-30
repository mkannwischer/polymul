#!/usr/bin/env python3
"""Random tests for all the algorithms."""
import subprocess

subprocess.run(["pydocstyle *.py"], shell=True, check=True, stderr=subprocess.PIPE)

tests = ["./02karatsuba.py", "./03toom.py", "./04ntt.py", "./05fft.py",
         "./06radix3fft.py", "./07incomplete.py", "./08goods.py"]

iterations = 1000

for test in tests:
    for _ in range(iterations):
        subprocess.check_output([test])
    print(f"ran {iterations} iterations of {test} test. ALL GOOD.")


