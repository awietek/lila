include options.mk
include sources.mk

rm     := rm -f
mkdir  := mkdir -p

objects = $(subst .cpp,.o,$(sources))
depends = $(subst .cpp,.d,$(sources))

depflags = -MT $@ -MMD -MP -MF $*.d


.PHONY: test 
test:  $(objects)
	$(cc) $(ccopt) $(ccarch) $(depflags) $(libraries) $(objects) $(testobjects) -o test/tests 

$(depends):
include $(depends)

# different building suites

.PHONY: clean
clean:
	$(RM) -r $(objects) $(mpiobjects) $(appobjects) $(testobjects) $(depends) $(appdepends) $(testdepends)

.PHONY: rebuild
rebuild: clean test

%.o: %.cpp %.d
%.o: %.cpp 
	$(cc) $(ccopt) $(ccarch) $(depflags) -c $< -o $@ $(includes)
