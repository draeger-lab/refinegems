FROM python:3.11-slim

WORKDIR /rg_cont

RUN apt-get update && apt-get install -y build-essential git && rm -rf /var/lib/apt/lists/*

ARG REFINEGEMS_EXTRAS="optional"
ARG INSTALL_EXTERNAL_TOOLS=true

COPY pyproject.toml ./
COPY README.md ./
COPY LICENSE THIRD_PARTY_LICENSES.md ./
COPY docs/source/images/logos/LOGO_LICENSE.md ./docs/source/images/logos/LOGO_LICENSE.md
COPY docs/source/images/GRAPHICS_LICENSE.md ./docs/source/images/GRAPHICS_LICENSE.md
COPY src/ ./src/

# Install the package using pip. Optional dependency groups can be adjusted with
# --build-arg REFINEGEMS_EXTRAS=chebi,ols,sbo or disabled with REFINEGEMS_EXTRAS=
RUN python -m pip install --no-cache-dir --upgrade pip \
    && if [ -n "$REFINEGEMS_EXTRAS" ]; then \
        python -m pip install --no-cache-dir ".[${REFINEGEMS_EXTRAS}]"; \
    else \
        python -m pip install --no-cache-dir .; \
    fi

# Install optional connected tools that are not distributed via the base package.
RUN if [ "$INSTALL_EXTERNAL_TOOLS" = "true" ]; then \
        python -m pip install --no-cache-dir "model-polisher@git+https://github.com/draeger-lab/MPClient" \
        && python -m pip install --no-cache-dir "masschargecuration@git+https://github.com/draeger-lab/MassChargeCuration" \
        && python -m pip install --no-cache-dir "bofdat@git+https://github.com/draeger-lab/BOFdat"; \
    fi

# Entrypoint for the container -> If set all CLI commands can be directly used
ENTRYPOINT ["refinegems"]
CMD ["-h"]
