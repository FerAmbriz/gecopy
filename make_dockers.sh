#!/bin/bash

folders=('cnvkit' 'exomedepth' 'panelcnmops')

for i in ${folders[@]};do
  mv $i/Dockerfile .
  docker build -t ambrizbiotech/$i .
  mv Dockerfile $i
done
