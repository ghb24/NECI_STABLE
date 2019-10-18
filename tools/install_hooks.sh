#!/usr/bin/env bash


PROJECT_SOURCE_DIR=$1

GIT_DIR=$(cd $PROJECT_SOURCE_DIR ; (cd $(git rev-parse --git-dir) ; pwd))

mkdir -p "$GIT_DIR/hooks"

DEST="$GIT_DIR/hooks/"

cp "$PROJECT_SOURCE_DIR/tools/pre-commit" $DEST
