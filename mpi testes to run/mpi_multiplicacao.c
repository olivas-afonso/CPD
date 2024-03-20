// Program 12.8 Release 1.3 on 30.08.2006
// The Parallel Gauss-Seidel Method (checkerboard decomposion)
  bool IsPrintEnabled;
  double *u;
  char methodName[] = "MPI Zeidel Block";
  TParams params;
  int mpiSsize;
  MPI_Datatype BlockType;
  MPI_Datatype VectorType;

int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);

  if (argc < 2)
  {
    printf("USAGE: GZ_MPIB.exe N\n");
    return 1;
  }
  if (argc == 2)
    IsPrintEnabled = 0;
  else
    IsPrintEnabled = atoi(argv[2]);

  MPI_Comm_rank(MPI_COMM_WORLD, ¶ms.mpiRank);
  MPI_Comm_size(MPI_COMM_WORLD, ¶ms.mpiSize);

  if (params.mpiRank)
    IsPrintEnabled = 0;

  if (!params.mpiRank)
    printf("Hello from %s method\n", methodName);

  int N = atoi(argv[1]);

  params.matrixWidth = N;
  params.matrixHeight = N;
  params.delta = 0.1;
  params.mainMatrix = CreateMatrix(N);

  GZ_Par(N, params);  

  MPI_Finalize();

  return 0;
}

int GZ_Par(int N, TParams params)
{
  if (initMethod(params))
  {
    printf("Error in initialization. Exit\n");
    return NULL;
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  double timeStart = MPI_Wtime();
  double currDelta, delta;
  
  // Main Loop
  int steps;
  for (steps = 0; steps < MAX_STEPS; steps++ )
  {
    currDelta = zeidelIteration(params);
    exchangeData(params);
    MPI_Allreduce(&currDelta, &delta, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    if (delta <= params.delta)
      break;
  }
  getData(params);
  MPI_Barrier(MPI_COMM_WORLD);
  double timeFinish = MPI_Wtime();

  if (!params.mpiRank && IsPrintEnabled)
    PrintMatrix(params.mainMatrix, N);

  deinitMethod(params);

  if (!params.mpiRank)
    printf("Results:\r\n- current eps = %f\r\n- required eps = %f\r\n- matrix Size = %d\r\n- iterations = %d\r\n- time = %f\r\n",
      delta, params.delta, N, steps, timeFinish - timeStart);

  return 0;
}

// init methods data, types etc.
int initMethod(TParams params)
{
  mpiSsize = 0;

  for (int i = 0; i < 10; i++)
    if (i * i + 1 == params.mpiSize)
      mpiSsize = i;
  if (!mpiSsize)
  {
    printf("Sorry, yours machine's count (%d != x*x+1) is not supports\n", params.mpiSize);
    return 1;
  }
  
  if (!params.mpiRank)
  {
    int BW = (params.matrixWidth - 2) / mpiSsize + 2;
    int BH = (params.matrixHeight - 2) / mpiSsize + 2;
    TParams pars;
    pars.matrixWidth = BW;
    pars.matrixHeight = BH;
    pars.mainMatrix = NULL;
    pars.anotherMatrix = NULL;
    pars.delta = params.delta;
    pars.mpiRank = -1;
    pars.mpiSize = -1;
    for(int i = 1; i < params.mpiSize; i++)
      MPI_Send(&pars, sizeof(pars), MPI_BYTE, i, 1, MPI_COMM_WORLD);
    
    MPI_Type_vector(BH, BW, params.matrixWidth, MPI_DOUBLE, &BlockType);
    MPI_Type_commit(&BlockType);
    int c = 1;
    for (i = 0; i < mpiSsize; i++)
      for (int j = 0; j < mpiSsize; j++)
        MPI_Send(params.mainMatrix + i * params.matrixWidth * (BH - 2) +
        j * (BW - 2), 1, BlockType, c++, 2, MPI_COMM_WORLD);
  }
  else
  {
    MPI_Status status;
    MPI_Recv(¶ms, sizeof(params), MPI_BYTE, 0, 1, MPI_COMM_WORLD, &status);
    params.mainMatrix = new double [params.matrixWidth * params.matrixHeight];
    MPI_Recv(params.mainMatrix, params.matrixWidth * params.matrixHeight,
      MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, &status);
    MPI_Type_vector(params.matrixHeight, 1, params.matrixWidth, MPI_DOUBLE, &VectorType);
    MPI_Type_commit(&VectorType);
  }
  return 0;
}

// exchange with edges of matrix
int exchangeData(TParams params)
{
  int c = 1;
  MPI_Status status;
  for (int i = 0; i < mpiSsize; i++)
    for (int j = 0; j < mpiSsize; j++)
    {
      if (params.mpiRank == c++)
      {
        int right  = (j == mpiSsize - 1)? MPI_PROC_NULL : i * mpiSsize + j + 1 + 1;
        int left   = (j == 0)? MPI_PROC_NULL : i * mpiSsize + j - 1 + 1;
        int bottom = (i == mpiSsize - 1)? MPI_PROC_NULL : (i + 1) * mpiSsize + j + 1;
        int top    = (i == 0)? MPI_PROC_NULL : (i - 1) * mpiSsize + j + 1;
        MPI_Sendrecv(params.mainMatrix + params.matrixWidth * (params.matrixHeight - 2),
          params.matrixWidth, MPI_DOUBLE, bottom, 4, params.mainMatrix, params.matrixWidth,
          MPI_DOUBLE, top, 4, MPI_COMM_WORLD, &status);
        MPI_Sendrecv(params.mainMatrix + params.matrixWidth, params.matrixWidth, MPI_DOUBLE, top, 5,
          params.mainMatrix + (params.matrixHeight - 1) * params.matrixWidth, params.matrixWidth,
          MPI_DOUBLE, bottom, 5, MPI_COMM_WORLD, &status);
        MPI_Sendrecv(params.mainMatrix + params.matrixWidth - 2, 1, VectorType, right, 6,
          params.mainMatrix, 1, VectorType, left, 6, MPI_COMM_WORLD, &status);
        MPI_Sendrecv(params.mainMatrix + 1, 1, VectorType, left, 7, params.mainMatrix +
          params.matrixWidth - 1, 1, VectorType, right, 7, MPI_COMM_WORLD, &status);
      }
    }

  return 0;
}

// Zeidel iteration
int zeidelIteration(TParams params)
{
  double * A = params.mainMatrix;

  //    printf("Seidel Iteration");
  if (params.mpiRank)
  {
    double OldValue, CurrentAccuracy, Accuracy = 0.0;
    for (int i = 1; i < params.matrixHeight - 1; i++)
      for (int j = 1; j < params.matrixWidth - 1; j++)
      {
        OldValue = A[i * params.matrixWidth + j];
        A[i * params.matrixWidth + j] = 0.25 * (A[(i - 1) * params.matrixWidth + j]
          + A[(i + 1) * params.matrixWidth + j] + A[i * params.matrixWidth + j - 1]
          + A[i * params.matrixWidth + j + 1]);
        CurrentAccuracy = fabs(A[i * params.matrixWidth + j] - OldValue);
        if (Accuracy < CurrentAccuracy)
          Accuracy = CurrentAccuracy;
      }
      return Accuracy;
  }
  return 0;
}

// delete matrix
int deinitMethod(TParams params)
{
  if (params.mpiRank)
  {
    MPI_Type_free(&VectorType);
    delete [] params.mainMatrix;
  }
  else
    MPI_Type_free(&BlockType);
  return 0;
}

// receive all computed data from slaves to master
int getData(TParams params)
{
  if (!params.mpiRank)
  {
    MPI_Status status;
    int BW = (params.matrixWidth - 2) / mpiSsize + 2;
    int BH = (params.matrixHeight - 2) / mpiSsize + 2;
    int c = 1;
    for (int i = 0; i < mpiSsize; i++)
      for (int j = 0; j < mpiSsize; j++)
        MPI_Recv(params.mainMatrix + i * params.matrixWidth * (BH - 2) +
        j * (BW - 2), 1, BlockType, c++, 1234, MPI_COMM_WORLD, &status);
  }
  else
  {
    MPI_Send(params.mainMatrix, params.matrixWidth * params.matrixHeight,
      MPI_DOUBLE, 0, 1234, MPI_COMM_WORLD);
  }

  return 0;
}

