#include <string.h>
#include <stdio.h>
#include <Omega_h_fail.hpp>
#include <egads.h>

#define CALL(f) OMEGA_H_CHECK(EGADS_SUCCESS == (f))

int main(int argc, char** argv) {
  if (argc != 3) {
    Omega_h_fail("usage: ./egads2lite model-in model-out\n");
  }

  size_t nbytes;
  ego context, model;
  char *stream;
  FILE *fp;

  CALL(EG_open(&context));
  CALL(EG_loadModel(context, 0, argv[1], &model));

  CALL(EG_exportModel(model, &nbytes, &stream));
  
  fp = fopen(argv[2], "wb");
  if (fp == NULL) exit(EXIT_FAILURE);
  fwrite(stream, sizeof(char), nbytes, fp);
  fclose(fp);
  EG_free(stream);

  CALL(EG_deleteObject(model));
  CALL(EG_close(context));
  return 0;
}
