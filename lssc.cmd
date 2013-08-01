SET pathfile=C:\Users\Yohann\Documents\GitHub\LSSC

SET execfile=%pathfile%\LSSC_build\bin\Release\Main.exe

SET name=lenaBW
SET extension=png
:: ARG 1 ==> path for original picture
SET img=%pathfile%\Images\%name%.%extension%
:: ARG 2 ==> noise standard deviation
SET sigma=25.0
:: ARG 3..7 ==> path for picture saving
SET imgNoisy=%pathfile%\Images\%name%_noisy.%extension%
SET imgFinal=%pathfile%\Images\%name%_final.%extension%
SET imgDiff=%pathfile%\Images\%name%_diff.%extension%
SET imgBias=%pathfile%\Images\%name%_bias.%extension%
SET imgDiffBias=%pathfile%\Images\%name%_diff_bias.%extension%
:: ARG 8 ==> boolean for bias
SET doBias=0

SET cmd=%execfile% %img% %sigma% %imgNoisy% %imgFinal% %imgDiff% %imgBias% %imgDiffBias% %doBias%

start %cmd%