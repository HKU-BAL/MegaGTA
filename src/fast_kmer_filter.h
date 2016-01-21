#include <stdlib.h>
#include <stdio.h>
#include <string.h>


class FastKmerFilter {
  private:
    struct RefKmer {
      private:
        int model_pos;
        int ref_file_index;
        string ref_seq_id;
      public:
        bool equals(RefKmer *k);
        int hashCode();
    };

  private:
    static void processSeq();

  public:
    static int find_start(int argc, char **argv);
};

