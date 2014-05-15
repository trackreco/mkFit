#ifdef ZMAJLA_ZAJLA

  std::cout << "LLLLL\n";

  for (int i = 0; i < M; ++i) { for (int j = 0; j < M; ++j)
      printf("%8f ", mul[0](i,j)); printf("\n");
  } printf("\n");
  for (int i = 0; i < M; ++i) { for (int j = 0; j < M; ++j)
      printf("%8f ", mpl.At(i,j,0)); printf("\n");
  } printf("\n");
  for (int i = 0; i < M; ++i) { for (int j = 0; j < M; ++j)
      printf("%8f ", muz[0](i,j)); printf("\n");
  } printf("\n");
  for (int i = 0; i < M; ++i) { for (int j = 0; j < M; ++j)
      printf("%8f ", mpz.At(i,j,0)); printf("\n");
  } printf("\n");

  for (int i = 0; i < N; ++i)
    for (int j = 0; j < M; ++j)
      for (int k = 0; k < M; ++k)
      {
        if (mul[i](j,k) != mpl.At(j, k, i))
          std::cout << i <<" "<< j <<" "<< k <<" "<< mul[i](j,k) <<" "<< mpl.At(j, k, i) << "\n";
      }

  std::cout << "ZZZZZ\n";

  for (int i = 0; i < N; ++i)
    for (int j = 0; j < M; ++j)
      for (int k = 0; k < M; ++k)
      {
        if (muz[i](j,k) != mpz.At(j, k, i))
          std::cout << i <<" "<< j <<" "<< k <<" "<< muz[i](j,k) <<" "<< mpz.At(j, k, i) << "\n";
      }

  std::cout << "RESRES\n";

#endif
