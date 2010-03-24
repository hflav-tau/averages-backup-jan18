SUBDIRS = $(subst /Makefile,,$(wildcard */Makefile))

.PHONY: subdirs $(SUBDIRS) combos pdg update_log plot

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

combos pdg update_log plot:
	for dir in $(SUBDIRS); do \
	  echo  $(MAKE) -C $$dir $@; \
	  $(MAKE) -C $$dir $@; \
	done
