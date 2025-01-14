import setuptools
import os
import eisd

here = os.path.abspath(os.path.dirname(__file__))

# Get the long description from the README file
# with open(path.join(here, 'README.md')) as f:
    # long_description = f.read()

if __name__ == "__main__":
    setuptools.setup(
        name='eisd',
        author='James Lincoff, Oufan Zhang, THGLab',
        project_urls={
            'Source': 'https://github.com/Oufan75/X-EISD'
        },
        description=
        'Experimental Inferential Structure Determination',
        keywords=[
            'Bayesian Optimization', 'Intrinsically Disordered Proteins'
        ],
        license='MIT',
        packages=setuptools.find_packages(),
        include_package_data=True,

        install_requires=[
            'numpy', 'pandas'
        ],
        extras_require={
            'docs': [
                'sphinx',
                'sphinxcontrib-napoleon',
                'sphinx_rtd_theme',
                'numpydoc',
                'nbsphinx'
            ],
            'tests': [
                'pytest',
                'pytest-cov',
                'pytest-pep8',
                'tox',
            ],
        },
        tests_require=[
            'pytest',
            'pytest-cov',
            'pytest-pep8',
            'tox',
        ],
        classifiers=[
            'Development Status :: 4 - Beta',
            'Natural Language :: English',
            'Intended Audience :: Science/Research',
            'Programming Language :: Python :: 3.7',
        ],
        zip_safe=False,
    )
