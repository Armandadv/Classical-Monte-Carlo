#!/usr/bin/env bash
# Build script (plain javac, no Maven/Gradle).
#   ./build.sh          -> compiles src/ into bin/ and test/ into bin_test/
#   ./build.sh src      -> compiles only src/
# Draw.java is excluded: it needs Processing 4 (circle()), which is not on Maven Central;
# the IDE (VSCode/Eclipse JDT) still compiles it when Processing 4 core.jar is on the class path.
set -euo pipefail
cd "$(dirname "$0")"
CP="lib/*"
# Separateur de classpath : ';' sur Windows (javac.exe), ':' ailleurs.
case "$(uname -s)" in CYGWIN*|MINGW*|MSYS*) SEP=";" ;; *) SEP=":" ;; esac
mkdir -p bin bin_test
SRC=$(ls src/*.java | grep -v '/Draw.java$')
javac -Xlint:-options -d bin -cp "$CP" $SRC
if [[ "${1:-all}" == "all" && -d test ]] && ls test/*.java >/dev/null 2>&1; then
  javac -Xlint:-options -d bin_test -cp "bin$SEP$CP" test/*.java
fi
echo "build OK"
