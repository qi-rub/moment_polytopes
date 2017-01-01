.PHONY: test docs upload-docs

test:
	SAGE_PATH=. sage -python -m pytest --doctest-glob="*.rst" --tb=short -s -vvvvv

docs:
	SAGE_PATH=. make -C docs html

upload-docs:
	git subtree push --prefix docs/_build/html origin gh-pages
