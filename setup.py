import setuptools

setuptools.setup(
    name="Koester_DAmodels",
    version="0.1.0",
    url="https://github.com/janvanroestel/Koester_DAmodels",
    author="Jan van Roestel",
    author_email="jcjvanroestel@gmail.com",
    description="Koestel DA white dwarf models with interpolator",
    long_description=open('README.md').read(),
    packages=setuptools.find_packages(),
    install_requires=['setuptools-git'],
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
    ],
    include_package_data=True,
    package_data={'': ['models_koester2_1571804487/*dat.txt']},
)
