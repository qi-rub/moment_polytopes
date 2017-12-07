import ast, re, io, os.path
from setuptools import setup

# determine version (adapted from mitsuhiko)
VERSION_RE = re.compile(r'__version__\s+=\s+(.*)')
with open('moment_polytopes/__init__.py', 'rb') as f:
    version = VERSION_RE.search(f.read().decode('utf-8')).group(1)
    version = str(ast.literal_eval(version))

# read long description
long_description = io.open(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), 'README'),
    encoding='utf-8').read()

setup(
    name='moment_polytopes',
    version=version,
    description='A SageMath package for computing moment polytopes.',
    long_description=long_description,
    author='Michael Walter',
    author_email='michael.walter@gmail.com',
    url='https://github.com/catch22/moment_polytopes',
    license='MIT',
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Intended Audience :: Science/Research'
        'Topic :: Scientific/Engineering :: Mathematics',
        'Development Status :: 3 - Alpha',
        'Programming Language :: Python :: 2',
    ],
    packages=['moment_polytopes'],
    install_requires=['tabulate'],
    extras_require={
        'dev': [
            'bibtex-pygments-lexer',
            'pytest>=3.3',
            'Sphinx',
        ]
    },
)
