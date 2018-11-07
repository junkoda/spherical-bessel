all: check

.PHONY: check clean coverage

check:
	python3 test.py

# pip3 install coverage
# `open` opens browser in Mac
coverage:
	coverage run test.py
	coverage html --omit='*numpy*'
	open htmlcov/index.html

clean:
	rm -r .coverage htmlcov
