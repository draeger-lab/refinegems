FROM python:3.14-slim

WORKDIR /rg_cont

RUN apt-get update && apt-get install -y build-essential git && rm -rf /var/lib/apt/lists/*

ARG INSTALL_EXTERNAL_TOOLS=true

COPY pyproject.toml ./
COPY README.md ./
COPY LICENSE THIRD_PARTY_LICENSES.md ./
COPY docs/source/images/logos/LOGO_LICENSE.md ./docs/source/images/logos/LOGO_LICENSE.md
COPY docs/source/images/GRAPHICS_LICENSE.md ./docs/source/images/GRAPHICS_LICENSE.md
COPY src/ ./src/

# Install the package with all runtime extras. The docs extra is intentionally
# excluded from the Docker image.
RUN python -m pip install --no-cache-dir --upgrade pip \
    && python -m pip install --no-cache-dir ".[optional]"

# Install optional connected tools that are not distributed via the base package.
RUN if [ "$INSTALL_EXTERNAL_TOOLS" = "true" ]; then \
        python -m pip install --no-cache-dir "model-polisher@git+https://github.com/draeger-lab/MPClient" \
        && python -m pip install --no-cache-dir "masschargecuration@git+https://github.com/draeger-lab/MassChargeCuration" \
        && python -m pip install --no-cache-dir "bofdat@git+https://github.com/draeger-lab/BOFdat"; \
    fi

# Entrypoint for the container -> If set all CLI commands can be directly used
ENTRYPOINT ["refinegems"]
CMD ["-h"]
