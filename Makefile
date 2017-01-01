.PHONY: test docs upload-docs

test:
docs:
	SAGE_PATH=. make -C docs html

upload-docs:
	git subtree push --prefix docs/_build/html origin gh-pages
