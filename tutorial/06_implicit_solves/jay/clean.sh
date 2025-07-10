#!/bin/bash

#user error
read -p "you sure you want to clean jay??: " choice

if [ "$choice" == "yes" ]; then
    rm -rf src/ *.txt bin/ __pycache__ .vscode main profile profile.dSYM
else
    echo "abort cleanup"
fi

# rm -rf src/ *.txt bin/ __pycache__

