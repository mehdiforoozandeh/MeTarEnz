# Extends the original 2020 image with procps, so workflow managers such as
# Nextflow can track processes with `ps`.
# The base is Debian buster, which is end-of-life: apt must use archive.debian.org.
FROM mforooz/metarenz:latest

RUN echo "deb http://archive.debian.org/debian buster main" > /etc/apt/sources.list \
 && apt-get -o Acquire::Check-Valid-Until=false update \
 && apt-get install -y --no-install-recommends procps \
 && rm -rf /var/lib/apt/lists/*
