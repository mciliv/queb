
export const isMobileDevice = () =>
  typeof window.matchMedia === 'function'
    ? window.matchMedia('(pointer: coarse) and (hover: none)').matches
    : false;
