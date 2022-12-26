ARG VERSION="0.4"

FROM alpine:3.15
ARG VERSION
RUN apk update && apk add --no-cache wget make gcc zlib-dev libc-dev mandoc man-pages
RUN wget -q https://github.com/lh3/miniprot/releases/download/v${VERSION}/miniprot-${VERSION}.tar.bz2 && \
    tar xjf miniprot-${VERSION}.tar.bz2 && \
    cd miniprot-${VERSION} && \
    make -j4

FROM alpine:3.15
ARG VERSION
RUN apk add mandoc
COPY --from=0 miniprot-${VERSION}/miniprot /usr/local/bin/miniprot
COPY --from=0 miniprot-${VERSION}/miniprot.1 /usr/share/man/man1/

WORKDIR /data
ENTRYPOINT ["miniprot"]

## Build:
# docker build --no-cache -t miniprot .
## Run: 
# docker run --rm miniprot
## Run as regular user mounting pwd: 
# docker run --rm -v $(pwd):$(pwd) -w $(pwd) -u $(id -u):$(id -g) miniprot
## Print man page
# docker run --rm --entrypoint man miniprot miniprot