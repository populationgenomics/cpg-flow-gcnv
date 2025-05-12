FROM australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_hail:0.2.134.cpg1 AS basic

ENV PYTHONDONTWRITEBYTECODE=1

# now do some fun stuff, installing ClinvArbitration
WORKDIR /cpg_flow_gcnv

COPY src src/
COPY pyproject.toml README.md ./

# pip install but don't retain the cache files
RUN pip install --no-cache-dir .
