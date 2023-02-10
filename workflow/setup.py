from setuptools import find_packages, setup


requirements = []


setup(
    name='absplice_scripts',
    packages=find_packages(),
    version='0.1.0',
    description='Scripts related with AbSplice ',
    author='Nils Wagner and M.Hasan Celik',
    license='MIT',
    test_suite='tests',
    tests_require=['pytest'],
    install_requires=requirements
)
