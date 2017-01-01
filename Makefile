.PHONY: upload-docs

test:

upload-docs:
	git subtree push --prefix docs/_build/html origin gh-pages
