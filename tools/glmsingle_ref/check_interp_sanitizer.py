#!/usr/bin/env python3
"""Regression check for the private GLMsingle interpolator's native bounds.

Run from any directory with a sanitizer-capable compiler, for example:
    CXX=clang++ python3 tools/glmsingle_ref/check_interp_sanitizer.py

Compile the function extracted verbatim from production with ASan and UBSan.
The R tests exercise the complete exported alpha kernel; this harness makes
invalid native reads observable even when their returned value looks like NaN.
"""

import os
from pathlib import Path
import re
import shlex
import subprocess
import tempfile


def main():
    root = Path(__file__).resolve().parents[2]
    source = (root / "src" / "glmsingle_kernels.cpp").read_text()
    match = re.search(r"^static double np_interp\(.*?^\}", source, re.M | re.S)
    if match is None:
        raise RuntimeError("Cannot locate the production np_interp function")
    harness = r"""
#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <vector>
#define NA_REAL std::numeric_limits<double>::quiet_NaN()
""" + match.group() + r"""
int main() {
    std::vector<double> bad(10, NA_REAL), fp(10, 1.0);
    assert(std::isnan(np_interp(0.5, bad.data(), fp.data(), 10)));
    assert(std::isnan(np_interp(0.5, nullptr, nullptr, 0)));
    double x[] = {0.0, 0.5, 1.0}, f[] = {0.0, 2.0, 4.0};
    assert(np_interp(0.25, x, f, 3) == 1.0);
    assert(np_interp(-1.0, x, f, 3) == 0.0);
    assert(np_interp(2.0, x, f, 3) == 4.0);
    double ties[] = {0.5, 0.5, 1.0};
    assert(np_interp(0.5, ties, f, 3) == 2.0);
    assert(np_interp(0.5, x, f, 1) == 0.0);
}
"""
    with tempfile.TemporaryDirectory(prefix="glms-interp-sanitizer-") as directory:
        path = Path(directory)
        src = path / "check.cpp"
        binary = path / "check"
        src.write_text(harness)
        compiler = shlex.split(os.environ.get("CXX", "clang++"))
        subprocess.run(compiler + ["-std=c++11", "-O1", "-g",
                                   "-fsanitize=address,undefined", str(src),
                                   "-o", str(binary)], check=True)
        subprocess.run([str(binary)], check=True)
    print("GLMsingle interpolation bounds: ASan/UBSan passed")


if __name__ == "__main__":
    main()
