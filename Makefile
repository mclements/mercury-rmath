
install: rmath.m
	rm -rf Mercury lib # is this required?
	mmc --make  --c-include-directory `pkg-config --cflags libRmath` `pkg-config --libs libRmath` --no-libgrade --libgrade asm_fast.gc --install-prefix . -E librmath.install

test: install test_rmath.m
	mmc --make --mld ./lib/mercury --ml rmath `pkg-config --libs libRmath` test_rmath && ./test_rmath

PROXY += test install
