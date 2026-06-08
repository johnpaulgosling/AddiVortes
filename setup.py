from __future__ import annotations

import sys

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext


class BuildExt(build_ext):
    c_opts = {
        "msvc": ["/std:c++20", "/O2"],
        "unix": ["-std=c++20", "-O3"],
    }

    def build_extensions(self) -> None:
        compiler_type = self.compiler.compiler_type
        opts = self.c_opts.get(compiler_type, [])
        if compiler_type == "unix" and sys.platform == "darwin":
            opts = [*opts, "-stdlib=libc++", "-mmacosx-version-min=10.15"]
        for ext in self.extensions:
            ext.extra_compile_args = [*getattr(ext, "extra_compile_args", []), *opts]
        super().build_extensions()


def ext_modules() -> list[Extension]:
    import pybind11

    return [
        Extension(
            "addivortes._core",
            ["python_src/addivortes_python.cpp"],
            include_dirs=[pybind11.get_include()],
            language="c++",
        )
    ]


setup(cmdclass={"build_ext": BuildExt}, ext_modules=ext_modules())
