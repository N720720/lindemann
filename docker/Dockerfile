FROM python:3.7-slim-buster

ENV LANG=C.UTF-8 \
  LC_ALL=C.UTF-8 \
  PATH="${PATH}:/root/.poetry/bin"

RUN apt-get update && \
  apt-get install -y --no-install-recommends \
  curl \
  make \
  && rm -rf /var/lib/apt/lists/*

COPY Makefile pyproject.toml ./

RUN make download-poetry && \
  make install NO_PRE_COMMIT=1 && \
  make build && \
  rm Makefile pyproject.toml

CMD mkdir -p /workspace

RUN cp /dist /workspace

WORKDIR /workspace

RUN  python -m pip install dist/lindemann-0.3.0-py3-none-any.whl
