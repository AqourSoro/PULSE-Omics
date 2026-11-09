from setuptools import setup, find_packages

setup(
    name="pulse",
    version="2.0.0",
    description="Patient Unified Longitudinal Signal Engine",
    author="PULSE Team",
    packages=find_packages(),
    install_requires=[
        "torch>=1.9.0",
        "numpy>=1.19.0",
        "tqdm>=4.60.0"
    ],
    python_requires=">=3.8",
)
