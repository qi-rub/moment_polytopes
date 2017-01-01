.PHONY: test docs upload-docs

TEST_FLAGS?=-s

test:
	SAGE_PATH=. sage -python -m pytest --tb=short $(TEST_FLAGS)

docs:
	SAGE_PATH=. make -C docs html

upload-docs:
	git subtree push --prefix docs/_build/html origin gh-pages
