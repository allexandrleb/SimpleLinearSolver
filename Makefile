
.SUFFIXES:
.SUFFIXES: .c .o

OBJDIR=tmp
ROOTDIR=$(pwd)

$(OBJDIR)/%.o: %.c | $(OBJDIR)
	gcc -c $< -o $@

all: $(OBJDIR)/nlf.o $(OBJDIR)/main.o $(OBJDIR)/linalg.o
	gcc $(OBJDIR)/nlf.o $(OBJDIR)/main.o $(OBJDIR)/linalg.o -o prog -lm -lgomp

$(OBJDIR):
	mkdir $(OBJDIR)

.PHONY : clean 
clean:
	rm prog
	rm $(OBJDIR)/*.o
	rm -r $(OBJDIR)
