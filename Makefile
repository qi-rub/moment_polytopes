.PHONY: test docs upload-docs

TEST_FLAGS?=-s

test:
	SAGE_PATH=. sage -python -m pytest --tb=short --doctest-modules moment_polytopes tests $(TEST_FLAGS)

docs:
	SAGE_PATH=. sage -ipython nbconvert examples/qmp.ipynb --to rst
	mv examples/qmp.rst docs/qmp.rst
	SAGE_PATH=. make -C docs html

upload-docs:
	git subtree push --prefix docs/_build/html origin gh-pages
	#git push origin `git subtree split --prefix docs/_build/html master`:gh-pages --force
