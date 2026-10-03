"use strict";

importScripts("inference.js");

self.onmessage = (event) => {
  try {
    const result = self.Inference.infer(event.data.spectrum, event.data.options);
    self.postMessage({ ok: true, result, requestId: event.data.requestId });
  } catch (error) {
    self.postMessage({ ok: false, error: error.message, requestId: event.data.requestId });
  }
};
