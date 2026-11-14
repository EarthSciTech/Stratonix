# Stratonix command-line program and C API in a container image.
#   docker run --rm -v "$PWD:/work" ghcr.io/earthscitech/stratonix:latest build model.yaml -o out
# Examples are in /opt/stratonix/share/stratonix/examples inside the image.
FROM ubuntu:24.04 AS build
RUN apt-get update \
 && apt-get install -y --no-install-recommends cmake ninja-build g++ pkg-config \
      libyaml-dev libcairo2-dev zlib1g-dev ca-certificates \
 && rm -rf /var/lib/apt/lists/*
COPY . /src
RUN cmake -S /src -B /build -G Ninja -DCMAKE_BUILD_TYPE=Release \
      -DSTRATONIX_BUILD_TESTS=OFF -DSTRATONIX_BUILD_BENCHMARKS=OFF \
 && cmake --build /build \
 && cmake --install /build --prefix /opt/stratonix

FROM ubuntu:24.04
RUN apt-get update \
 && apt-get install -y --no-install-recommends libyaml-0-2 libcairo2 zlib1g fonts-dejavu-core \
 && rm -rf /var/lib/apt/lists/*
COPY --from=build /opt/stratonix /opt/stratonix
ENV PATH=/opt/stratonix/bin:$PATH \
    LD_LIBRARY_PATH=/opt/stratonix/lib
LABEL org.opencontainers.image.title="Stratonix" \
      org.opencontainers.image.description="C++ toolkit for 2D/3D structural geomodelling: command-line program and C API" \
      org.opencontainers.image.source="https://github.com/EarthSciTech/Stratonix" \
      org.opencontainers.image.licenses="MIT"
WORKDIR /work
ENTRYPOINT ["stratonix"]
CMD ["version"]
