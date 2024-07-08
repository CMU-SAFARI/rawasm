
# This simple Makefile is used to install the rawasm extension for miniasm.
# It can either download miniasm and patch it, or patch an existing miniasm dir.
# This patch only applies on src and not binaries
# Librarie are pre-compiled but it is possible to build them again with the lib_rebuild rule

MINIASM_PATH = ./miniasm

# Download miniasm and install rawasm extension
all: get_miniasm install 

# Download vanilla miniasm
get_miniasm:
	git clone https://github.com/lh3/miniasm.git $(MINIASM_PATH)

# Install rawasm
install: patch build

# Patch miniasm with rawasm additional files
patch:
	mkdir $(MINIASM_PATH)/rawasm 
	cp -r include lib src $(MINIASM_PATH)/rawasm/
	cp patch/* $(MINIASM_PATH)/

# Compile miniasm with the new rawasm files
build:
	cd $(MINIASM_PATH); make

# Build the libraries (they come precompiled by default)
lib_rebuild:

# Remove the rawasm installation (miniasm is deleted as well)
clean:
	yes | rm -r $(MINIASM_PATH)

.PHONY: all patch clean
