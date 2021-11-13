#!/bin/fish
set -gx PYTHONPATH $PYTHONPATH (pwd)/core
echo "now python is aware of the path: "$PYTHONPATH
