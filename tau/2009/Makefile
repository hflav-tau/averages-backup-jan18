SUBDIRS = $(subst /Makefile,,$(wildcard */Makefile))

.PHONY: subdirs $(SUBDIRS) combos pdg update_log plot

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

combos pdg plot update_log update_pdglog update_plotlog update_all:
	for dir in $(SUBDIRS); do \
	  echo  $(MAKE) -C $$dir $@; \
	  $(MAKE) -C $$dir $@; \
	done
