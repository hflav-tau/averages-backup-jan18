SUBDIRS_DASH = $(subst /Makefile,,$(wildcard *-*/Makefile))
SUBDIRS = $(filter-out $(SUBDIRS_DASH),$(subst /Makefile,,$(wildcard */Makefile)))
TARGETS = all combos alucomb pdg plot update_log update_alucomblog update_pdglog update_plotlog update_all

.PHONY: info subdirs $(SUBDIRS) $(TARGETS)

info:
	@echo "please use targets: $(TARGETS)"
	@echo "sub-directories: $(SUBDIRS)"

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

$(TARGETS):
	for dir in $(SUBDIRS); do \
	  echo  $(MAKE) -C $$dir $@; \
	  $(MAKE) -C $$dir $@; \
	done
