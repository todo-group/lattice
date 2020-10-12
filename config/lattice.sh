#!/bin/sh

LATTICE_VERSION=0.1

rm -rf $BASEDIR/lattice
curl -L https://github.com/todo-group/lattice/archive/${LATTICE_VERSION}.tar.gz | tar zxf - -C ${BASEDIR}
mv -f ${BASEDIR}/lattice-${LATTICE_VERSION}/lattice ${BASEDIR}
mv -f ${BASEDIR}/lattice-${LATTICE_VERSION}/config/lattice.sh ${BASEDIR}/config
rm -rf ${BASEDIR}/lattice-${LATTICE_VERSION}
