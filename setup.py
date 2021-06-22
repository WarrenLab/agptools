"""Setup the package"""
import setuptools

setuptools.setup(
    name="agptools",
    packages=setuptools.find_packages(),
    scripts=["agptools", "break_contigs"],
    install_requires=["screed"],
)
