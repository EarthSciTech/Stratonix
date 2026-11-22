# The wheel carries the platform's C library (stratonix/_lib), so it is a platform wheel that
# serves every Python 3 version (py3-none-<platform>): no extension module is compiled.
from setuptools import setup
from setuptools.dist import Distribution

try:
    from wheel.bdist_wheel import bdist_wheel as _bdist_wheel

    class bdist_wheel(_bdist_wheel):
        def finalize_options(self):
            super().finalize_options()
            self.root_is_pure = False

        def get_tag(self):
            _, _, plat = super().get_tag()
            return "py3", "none", plat

    cmdclass = {"bdist_wheel": bdist_wheel}
except ImportError:
    cmdclass = {}


class BinaryDistribution(Distribution):
    def has_ext_modules(self):
        return True


setup(cmdclass=cmdclass, distclass=BinaryDistribution)
