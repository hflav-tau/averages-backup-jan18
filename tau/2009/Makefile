SUBDIRS_DASH = $(filter-out *-*,$(subst /Makefile,,$(wildcard *-*/Makefile)))
SUBDIRS = $(filter-out $(SUBDIRS_DASH),$(subst /Makefile,,$(wildcard */Makefile)))
TARGETS = combos pdg plot update_log update_pdglog update_plotlog update_all

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
