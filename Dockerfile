# syntax=docker/dockerfile:1
FROM debian:bookworm-slim AS rtl-driver

ARG RTL_SDR_COMMIT=aed0ea19f3a273370a13c9009b96313c75d54c7b
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential ca-certificates cmake git libusb-1.0-0-dev pkg-config \
    && rm -rf /var/lib/apt/lists/*
RUN git init /src \
    && git -C /src fetch --depth 1 https://github.com/rtlsdrblog/rtl-sdr-blog.git "${RTL_SDR_COMMIT}" \
    && git -C /src checkout --detach FETCH_HEAD \
    && cmake -S /src -B /build \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -DCMAKE_INSTALL_LIBDIR=lib \
        -DDETACH_KERNEL_DRIVER=ON \
        -DINSTALL_UDEV_RULES=OFF \
        -DENABLE_ZEROCOPY=OFF \
    && cmake --build /build --parallel 2 \
    && DESTDIR=/install cmake --install /build \
    && mkdir -p /install/usr/local/share/doc/rtl-sdr \
    && cp /src/COPYING /install/usr/local/share/doc/rtl-sdr/COPYING \
    && git -C /src archive --format=tar.gz HEAD > /install/usr/local/share/doc/rtl-sdr/source.tar.gz \
    && printf '%s\n' "${RTL_SDR_COMMIT}" > /install/usr/local/share/doc/rtl-sdr/source-commit

FROM debian:bookworm-slim AS runtime
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates libusb-1.0-0 python3 python3-numpy \
    && rm -rf /var/lib/apt/lists/*
COPY --from=rtl-driver /install/usr/local/ /usr/local/
RUN ldconfig
ENV PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
WORKDIR /app
COPY dab_ppm.py /app/dab_ppm.py
WORKDIR /data
ENTRYPOINT ["python3", "/app/dab_ppm.py"]
CMD ["--help"]
