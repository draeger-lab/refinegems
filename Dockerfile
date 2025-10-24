FROM python:3.11-slim

WORKDIR /rg_cont

RUN apt-get update && apt-get install -y build-essential git && rm -rf /var/lib/apt/lists/*

COPY pyproject.toml ./
COPY src/ ./src/
COPY README.md ./

# Install the package using pip (setuptools will be used automatically)
RUN pip install .

# Install dependencies that cannot be specified in pyproject.toml
RUN pip install "model-polisher@git+https://github.com/draeger-lab/MPClient"
RUN pip install "masschargecuration@git+https://github.com/draeger-lab/MassChargeCuration"
RUN pip install "bofdat@git+https://github.com/draeger-lab/BOFdat"

# Entrypoint for the container -> If set all CLI commands can be directly used
ENTRYPOINT ["refinegems"]
CMD ["-h"]