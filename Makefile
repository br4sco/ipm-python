.PHONY: coherent

MODULES = ipm_furuta

coherent:
	$(eval VERSION := $(shell grep '^__version__ = ".*"$$' setup.py))
	@if [ -z '$(VERSION)' ]; then \
	  echo "Missing version. Check the setup.py file."; \
	  exit 1; \
	 fi
	echo '$(VERSION)'

	sed -i 's|^__version__ = ".*"$$|$(VERSION)|g' $(foreach mod, $(MODULES), $(mod)/__init__.py)
