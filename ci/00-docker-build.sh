#!/bin/sh

set -e

set -o pipefail

# This wrapper is run on the hosts that build the containers. The goal is
# to minimize the context, so that silly little changes to the ci/ tree
# don't necessarily trigger a full rebuild of the containers.

IMAGE="$1"

tmp=$(mktemp -d /tmp/XXXXXXXXXX)
trap "rm -rf $tmp" EXIT

(cd $(dirname $0)/ ; rsync -a --files-from=- ./ "$tmp"/) <<EOF
000-functions.sh
001-environment.sh
00-prepare-docker.sh
utilities/ncpus.sh
EOF

ci/00-dockerfile.sh > "$tmp/Dockerfile"

docker build -t $IMAGE --cache-from $IMAGE:latest $tmp
