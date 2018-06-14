#!/usr/bin/env bash

brew cask uninstall oclint
brew tap oclint/formulae
brew install oclint gcc llvm
echo 'export PATH="/usr/local/opt/llvm/bin:$PATH"' >> ~/.bash_profile