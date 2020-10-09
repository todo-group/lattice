#!/bin/sh

BASEDIR=$(cd "$(dirname $0)"; pwd)
SCRIPTS=$(find ${BASEDIR}/setup -name "*.sh")
echo ${SCRIPTS}
for i in ${SCRIPTS}; do
  env BASEDIR=${BASEDIR} sh $i
done
