SUBDIRS_EXCLUDE = $(subst /00-exclude-flag.txt,,$(wildcard */00-exclude-flag.txt))
SUBDIRS = $(filter-out $(SUBDIRS_EXCLUDE),$(subst /Makefile,,$(wildcard */Makefile)))
TARGETS = all combos alucomb pdg plot update_log update_alucomblog update_pdglog update_plotlog update_all makelink

.PHONY: info subdirs $(SUBDIRS) $(TARGETS)

info:
	@echo "please use targets: $(TARGETS)"
	@echo "sub-directories: $(SUBDIRS)"

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

$(TARGETS):
	@echo "will make in subdirs: $(SUBDIRS)"
	@for dir in $(SUBDIRS); do \
	  echo  $(MAKE) -C $$dir $@; \
	  $(MAKE) -C $$dir $@; \
	done
	@echo "make done in subdirs: $(SUBDIRS)"
