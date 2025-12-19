# SPDX-License-Identifier: LGPL-3.0-or-later
"""
Setup script for genepy package.

This setup.py is required to force platform-specific wheel creation
when bundling pre-built binaries (libpython_interface.so, atdyn).
"""

from setuptools import setup

try:
    from wheel.bdist_wheel import bdist_wheel as _bdist_wheel

    class bdist_wheel(_bdist_wheel):
        """Custom bdist_wheel command that creates platform-specific wheels."""

        def finalize_options(self):
            _bdist_wheel.finalize_options(self)
            # Force platform wheel (not pure Python)
            self.root_is_pure = False

        def get_tag(self):
            python, abi, plat = _bdist_wheel.get_tag(self)
            # Use generic Python tag since we don't have compiled extensions
            # but we do have platform-specific binaries
            python = "py3"
            abi = "none"
            return python, abi, plat

except ImportError:
    bdist_wheel = None

cmdclass = {}
if bdist_wheel is not None:
    cmdclass["bdist_wheel"] = bdist_wheel

setup(cmdclass=cmdclass)
